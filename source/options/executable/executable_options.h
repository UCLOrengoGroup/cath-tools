/// \file
/// \brief The executable_options class header

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

#ifndef EXECUTABLE_OPTIONS_H_INCLUDED
#define EXECUTABLE_OPTIONS_H_INCLUDED

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include "common/type_aliases.h"

#include <string>
#include <vector>

namespace cath { namespace opts { class options_block; } }

namespace cath {
	namespace opts {

		/// \brief Provide ABC-interface for an executable's options handling
		///
		/// Concrete derived classes should add their options blocks to this via add_options_block().
		/// In practice, they should probably contain those options blocks and perform any super::add_options_block()
		/// calls in their ctors.
		///
		/// The options blocks must still exist when parse_options() is called.
		///
		/// Several methods (eg get_variables_map() and get_error_or_help_string()) can only be called
		/// once parse_options() has been called.
		///
		/// At present, this class is kept simpler by being one-use only, meaning that parse_options()
		/// can only be called once and add_options_block() cannot be called after that.
		///
		/// Current (descending) order of precedence:
		///  - Command line options
		///  - Global environment variables (eg var  : "CATH_TOOLS_PDB_PATH" )
		///  - Global configuration file    (   file : "cath-tools.conf"     )
		class executable_options {
		private:
			/// \brief This affects the way Boost Program Options formats the outputting of the options description
			///
			/// Its behaviour isn't 100% clear but setting this value to roughly the character-width of a
			/// modern terminal prevents it prematurely wrapping and making a very long, narrow output.
			static constexpr size_t                  DEFAULT_PROG_OPS_LINE_LENGTH = 200;

			static const     std::string             CATH_TOOLS_ENVIRONMENT_VARIABLE_PREFIX;
			static const     boost::filesystem::path CATH_TOOLS_CONF_FILE;
			static const     path_vec                CATH_TOOLS_CONF_FILE_SEARCH_PATH;

			/// \brief A list of pointers to the options blocks to be processed during parsing
			std::vector<options_block *> all_options_blocks;

			/// \brief A list of pointers to the options blocks to be processed during parsing
			std::vector<options_block *> visible_options_blocks;

			/// \brief Whether or not options have been parsed
			///
			/// This is used to prevent post-parsing methods being called before parse_options() is called
			/// and to prevent pre-parsing methods being called before
			bool processed_options = false;

			/// \brief A string that is populated with any error/help messages that arise
			///        during parsing. May be queried with get_error_or_help_string()
			std::string error_or_help_string;

			/// \brief The Boost program_options variable map which stores details of the parsing.
			///        This can be queried by the protected method get_variables_map() for checking
			///        whether options were defaulted(), etc.
			boost::program_options::variables_map vm;

			/// \brief Provide a name for the executable
			///
			/// This is used in error/help strings
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			virtual std::string do_get_program_name() const = 0;

			/// \brief Review all specified options and return a string containing any errors or a help string
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			virtual std::string do_update_error_or_help_string(const boost::program_options::options_description &) const = 0;

			virtual boost::program_options::positional_options_description get_positional_options();

		protected:
			std::string get_standard_usage_error_string() const;
			std::string get_program_name() const;
			const boost::program_options::variables_map & get_variables_map() const;

			void add_options_block(options_block &);

		public:
			virtual ~executable_options() noexcept = default;

			void parse_options(const int &,
			                   const char * const []);

			std::string get_error_or_help_string() const;
		};

	}
}

#endif
