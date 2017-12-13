/// \file
/// \brief The bootstrap_mode header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_SRC_TEST_TEST_PREDICATE_BOOTSTRAP_MODE_HPP
#define _CATH_TOOLS_SOURCE_SRC_TEST_TEST_PREDICATE_BOOTSTRAP_MODE_HPP

#include <string>

namespace cath {
	namespace test {

		/// \brief Whether a test file should be bootstrapped - ie replaced with the "got" content if mismatching
		enum class bootstrap_mode : char {
			ALWAYS, ///< Always bootstrap the files
			IF_ENV, ///< Bootstrap the files if the get_bootstrap_env_var() environment variable is set
			NEVER   ///< Never bootstrap the file
		};

		std::string get_bootstrap_env_var();

		bool bootstrap_env_var_is_set();

		bool should_overwrite(const bootstrap_mode &);

	} // namespace test
} // namespace cath

#endif
