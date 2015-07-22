/// \file
/// \brief The misc_help_version_options_block class header

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

#ifndef MISC_HELP_VERSION_OPTIONS_BLOCK_H_INCLUDED
#define MISC_HELP_VERSION_OPTIONS_BLOCK_H_INCLUDED

#include "options/options_block/options_block.h"

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

			/// \brief Whether help has been requested
			bool help    = false;

			/// \brief Whether version information has been requested
			bool version = false;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual opt_str do_invalid_string() const override final;

		public:
			virtual ~misc_help_version_options_block() noexcept = default;

			bool get_help() const;
			bool get_version() const;

			std::string get_help_string(const boost::program_options::options_description &,
			                            const std::string &,
			                            const std::string &) const;
			std::string get_version_string(const std::string &,
			                               const std::string &) const;

			static const std::string CATH_BINARIES_VERSION;
			static const std::string CATH_BINARIES_VERSION_DATE;

			static const std::string PO_HELP;
			static const std::string PO_VERSION;
		};
	}
}

#endif