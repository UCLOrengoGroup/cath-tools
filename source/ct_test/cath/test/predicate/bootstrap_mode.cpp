/// \file
/// \brief The bootstrap_mode class definitions

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

#include "bootstrap_mode.hpp"

#include <cstdlib>

using ::std::string;

/// \brief Get the environment variable used for activate test bootstrapping
string cath::test::get_bootstrap_env_var() {
	return "BOOTSTRAP_TESTS";
}

/// \brief Return whether the bootstrap environment variable is set
bool cath::test::bootstrap_env_var_is_set() {
	return ( getenv( get_bootstrap_env_var().c_str() ) != nullptr );
}

/// \brief Return whether a file should be written under the specified bootstrap_mode
bool cath::test::should_overwrite(const bootstrap_mode &prm_bootstrap_mode ///< The bootstrap_mode for the file
                                  ) {
	return (
		( prm_bootstrap_mode == bootstrap_mode::ALWAYS )
		||
		(
			( prm_bootstrap_mode == bootstrap_mode::IF_ENV )
			&&
			bootstrap_env_var_is_set()
		)
	);
}
