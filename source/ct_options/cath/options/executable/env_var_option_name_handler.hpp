/// \file
/// \brief The env_var_option_name_handler class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_ENV_VAR_OPTION_NAME_HANDLER_HPP
#define _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_ENV_VAR_OPTION_NAME_HANDLER_HPP

#include <boost/program_options.hpp>

#include <string>

namespace cath {
	namespace opts {

		/// \brief Functor to be passed to Boost Program Options' parse_environment() function to do
		///        the work of translating environment variable names to options names
		///
		/// This achieves two things that can't be achieved by just passing a prefix string
		/// to parse_environment (at least, not at Boost 1.54.0) :
		///  -# Converting underscores to hyphens
		///  -# Optionally allowing unrecognised options
		///     (by searching for each candidate option name and returning "" to parse_environment()
		///      for each that's unrecognised)
		///
		/// Example usage: \snippet executable_options.cpp Using env_var_option_name_handler
		class env_var_option_name_handler final {
			std::string prefix;
			bool allow_unknown;
			const boost::program_options::options_description &the_options;

		public:
			env_var_option_name_handler(std::string,
			                            const bool &,
			                            const boost::program_options::options_description & = boost::program_options::options_description());

			std::string operator()(const std::string &) const;
		};

		std::string option_of_environment_variable_and_prefix(const std::string &,
		                                                      const std::string &);

		std::string environment_variable_prefix_of_program_name(const std::string &);
		
	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_EXECUTABLE_ENV_VAR_OPTION_NAME_HANDLER_HPP
