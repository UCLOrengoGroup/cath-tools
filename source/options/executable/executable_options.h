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
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include "common/boost_addenda/program_options/set_opt_str_from_prog_opts_try.h"
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
			/// \brief The line legnth to use when rendering program options
			///
			/// Its behaviour isn't 100% clear but setting this value to roughly the character-width of a
			/// modern terminal prevents it prematurely wrapping and making a very long, narrow output.
			static constexpr size_t                  DEFAULT_PROG_OPS_LINE_LENGTH = 200;

			static const     std::string             CATH_TOOLS_ENVIRONMENT_VARIABLE_PREFIX;
			static const     boost::filesystem::path CATH_TOOLS_CONF_FILE;
			static const     path_vec                CATH_TOOLS_CONF_FILE_SEARCH_PATH;

			/// \brief A list of (references to) the options blocks to be processed during parsing
			std::vector<std::reference_wrapper<options_block>> all_options_blocks;

			/// \brief Whether or not options have been parsed
			///
			/// This is used to prevent post-parsing methods being called before parse_options() is called
			/// and to prevent pre-parsing methods being called before
			bool processed_options = false;

			/// \brief A string that is populated with any error/help messages that arise
			///        during parsing. May be queried with get_error_or_help_string()
			opt_str error_or_help_string;

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
			virtual opt_str do_get_error_or_help_string() const = 0;

			virtual boost::program_options::positional_options_description get_positional_options();

			virtual std::string do_get_help_prefix_string() const = 0;
			virtual std::string do_get_help_suffix_string() const = 0;
			virtual std::string do_get_overview_string() const = 0;

			template <typename FN>
			void prog_opts_try(opt_str &,
			                   FN &&,
			                   const opt_str & = boost::none);

		protected:
			std::string get_standard_usage_error_string() const;
			std::string get_program_name() const;
			std::string get_help_prefix_string() const;
			std::string get_help_suffix_string() const;
			std::string get_overview_string() const;
			const boost::program_options::variables_map & get_variables_map() const;

			void add_options_block(options_block &);

			static void add_all_options_to_description(boost::program_options::options_description &,
			                                           options_block &,
			                                           const size_t &);
			static void add_visble_options_to_description(boost::program_options::options_description &,
			                                              options_block &,
			                                              const size_t &);

		public:
			virtual ~executable_options() noexcept = default;

			void parse_options(const int &,
			                   const char * const []);

			const opt_str & get_error_or_help_string() const;
		};

		/// \brief Try a program options action and handle any exceptions that are thrown
		template <typename Func>
		void executable_options::prog_opts_try(opt_str        &arg_error_string, ///< The optional error string to update with a description of any errors that occur
		                                       Func          &&arg_function,     ///< The function to perform
		                                       const opt_str  &arg_parsing_phase ///< The phase in which this parsing is occurring (or none)
		                                       ) {
			common::set_opt_str_from_prog_opts_try(
				arg_error_string,
				std::forward<Func>( arg_function ),
				get_program_name() + ": " + ( arg_parsing_phase ? *arg_parsing_phase + " " : std::string{} ),
				"\n" + get_standard_usage_error_string()
			);
		}

		/// \brief Return a new instance of the specified type of executable_options with the specified options parsed into it
		template <typename T>
		T make_and_parse_options(const int          &argc,  ///< The main()-style argc parameter
		                         const char * const  argv[] ///< The main()-style argv parameter
		                         ) {
			static_assert(
				std::is_base_of<executable_options, T>::value,
				"make_and_parse_options() can only be used to make types that inherit from executable_options."
			);
			T new_options;
			new_options.parse_options( argc, argv );
			return new_options;
		}
	}
}

#endif
