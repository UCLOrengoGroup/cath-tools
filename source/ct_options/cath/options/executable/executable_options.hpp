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

#ifndef _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_EXECUTABLE_OPTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_EXECUTABLE_OPTIONS_HPP

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "cath/common/argc_argv_faker.hpp"
#include "cath/common/boost_addenda/program_options/set_opt_str_from_prog_opts_try.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/options/executable/parse_sources.hpp"
#include "cath/options/options_block/string_options_block.hpp"

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
			/// \brief The line length to use when rendering program options
			///
			/// Its behaviour isn't 100% clear but setting this value to roughly the character-width of a
			/// modern terminal prevents it prematurely wrapping and making a very long, narrow output.
			static constexpr size_t                  DEFAULT_PROG_OPS_LINE_LENGTH = 200;

			static const     std::string             CATH_TOOLS_ENVIRONMENT_VARIABLE_PREFIX;
			static const     ::std::filesystem::path CATH_TOOLS_CONF_FILE;

			static path_vec CATH_TOOLS_CONF_FILE_SEARCH_PATH();

			/// \brief A list of string_options_block objects that just serve to put strings
			///        the options usage message
			///
			/// They are stored here:
			///  * because all_options_blocks stores reference and
			///  * so that the clients can add strings without having to worry about string_options_block lifetimes
			///
			/// Using deque not vector because the all_options_blocks references mustn't get invalidated
			std::deque<string_options_block> string_obj_blocks;

			/// \brief A list of (references to) the options blocks to be processed during parsing
			std::vector<std::reference_wrapper<options_block>> all_options_blocks;

			/// \brief Whether options have been parsed
			///
			/// This is used to prevent post-parsing methods being called before parse_options() is called
			/// and to prevent pre-parsing methods being called before
			bool processed_options = false;

			/// \brief A string that is populated with any error/help messages that arise
			///        during parsing. May be queried with get_error_or_help_string()
			str_opt error_or_help_string;

			/// \brief The Boost program_options variable map which stores details of the parsing.
			///        This can be queried by the protected method get_variables_map() for checking
			///        whether options were defaulted(), etc.
			boost::program_options::variables_map vm;

			/// \brief Provide a name for the executable
			///
			/// This is used in error/help strings
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			[[nodiscard]] virtual std::string do_get_program_name() const = 0;

			/// \brief Review all specified options and return a string containing any errors or a help string
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			[[nodiscard]] virtual str_opt do_get_error_or_help_string() const = 0;

			virtual boost::program_options::positional_options_description get_positional_options();

			[[nodiscard]] virtual std::string do_get_help_prefix_string() const = 0;
			[[nodiscard]] virtual std::string do_get_help_suffix_string() const = 0;
			[[nodiscard]] virtual std::string do_get_overview_string() const    = 0;

			template <typename FN>
			void prog_opts_try( str_opt &, FN &&, const str_opt & = ::std::nullopt );

		  protected:
			[[nodiscard]] std::string                                  get_standard_usage_error_string() const;
			[[nodiscard]] std::string                                  get_program_name() const;
			[[nodiscard]] std::string                                  get_help_prefix_string() const;
			[[nodiscard]] std::string                                  get_help_suffix_string() const;
			[[nodiscard]] std::string                                  get_overview_string() const;
			[[nodiscard]] const boost::program_options::variables_map &get_variables_map() const;

			void add_options_block( options_block & );
			void add_string( std::string );

			static void add_all_options_to_description( boost::program_options::options_description &,
			                                            options_block &,
			                                            const size_t & );
			static void add_visble_options_to_description( boost::program_options::options_description &,
			                                               options_block &,
			                                               const size_t & );

		  public:
			executable_options() = default;
			virtual ~executable_options() noexcept = default;

			executable_options(const executable_options &) = default;
			/// \TODO: Come consistently more recent Clangs than the 3.9.0 used by Travis-CI,
			///        make this noexcept
			executable_options(executable_options &&) = default;
			executable_options & operator=(const executable_options &) = default;
			/// \TODO: Come consistently more recent Clangs than the 3.9.0 used by Travis-CI,
			///        make this noexcept
			executable_options & operator=(executable_options &&) = default;

			void parse_options(const int &,
			                   const char * const [],
			                   const parse_sources & = parse_sources::CMND_ENV_AND_FILE);

			[[nodiscard]] str_opt get_error_or_help_string() const;
		};

		/// \brief Try a program options action and handle any exceptions that are thrown
		template <typename Func>
		void executable_options::prog_opts_try(str_opt        &prm_error_string, ///< The optional error string to update with a description of any errors that occur
		                                       Func          &&prm_function,     ///< The function to perform
		                                       const str_opt  &prm_parsing_phase ///< The phase in which this parsing is occurring (or nullopt)
		                                       ) {
			common::set_opt_str_from_prog_opts_try(
				prm_error_string,
				std::forward<Func>( prm_function ),
				get_program_name() + ": " + ( prm_parsing_phase ? *prm_parsing_phase + " " : std::string{} ),
				"\n" + get_standard_usage_error_string()
			);
		}

		/// \brief Return a new instance of the specified type of executable_options with the specified options parsed into it
		template <typename T>
		T make_and_parse_options(const int           &argc,                                                ///< The main()-style argc parameter
		                         const char * const   argv[],                                              ///< The main()-style argv parameter
		                         const parse_sources &prm_parse_sources = parse_sources::CMND_ENV_AND_FILE ///< The sources from which options should be parsed
		                         ) {
			static_assert(
				std::is_base_of_v<executable_options, T>,
				"make_and_parse_options() can only be used to make types that inherit from executable_options."
			);
			T new_options;
			new_options.parse_options( argc, argv, prm_parse_sources );
			return new_options;
		}

		/// \brief Return a new instance of the specified type of executable_options with the specified options parsed into it
		template <typename T>
		T make_and_parse_options(const str_vec       &args,                                                ///< The arguments
		                         const parse_sources &prm_parse_sources = parse_sources::CMND_ENV_AND_FILE ///< The sources from which options should be parsed
		                         ) {
			argc_argv_faker my_argc_argv_faker( args );
			return make_and_parse_options<T>(
				my_argc_argv_faker.get_argc(),
				my_argc_argv_faker.get_argv(),
				prm_parse_sources
			);
		}
	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_EXECUTABLE_OPTIONS_HPP
